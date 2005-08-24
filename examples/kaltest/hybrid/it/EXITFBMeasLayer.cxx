//*************************************************************************
//* ===================
//*  EXITFBMeasLayer Class
//* ===================
//*
//* (Description)
//*   Sample measurement layer class used by EXITFBHit.
//* (Requires)
//* (Provides)
//*     class EXITFBMeasLayer
//* (Update Recored)
//*   2003/09/30  Y.Nakashima       Original version.
//*   2005/07/25  Kim, Youngim      Forward & Backward versoin.
//*
//*************************************************************************
//

#include "EXITFBMeasLayer.h"
#include "EXITFBHit.h"
#include "EXITKalDetector.h"
#include "TVTrack.h"
#include "TRandom.h"
#include <iostream>

using namespace std;
ClassImp(EXITFBMeasLayer)
                                                                                
EXITFBMeasLayer::EXITFBMeasLayer(TMaterial &min,
                                 TMaterial &mout,
                                 TVector3  &xc,
                                 Double_t   rin,
                                 Double_t   rout,
                                 Double_t   sigmax,
                                 Double_t   sigmay,
                                 Bool_t     type)
               : EXVMeasLayer(min, mout, type),
                 TPlane(xc),
                 fRin(rin),
                 fRout(rout),
                 fSigmaX(sigmax),
                 fSigmaY(sigmay)
{
}

EXITFBMeasLayer::~EXITFBMeasLayer()
{
}

TKalMatrix EXITFBMeasLayer::XvToMv(const TVector3 &xv) const
{
   // Calculate hit coordinate information:
   //	mv(0,0) = x
   //     (1,0) = y

   TKalMatrix mv(kMdim,1);
   mv(0,0)  = xv.X();
   mv(1,0)  = xv.Y();
   return mv;
}

TKalMatrix EXITFBMeasLayer::XvToMv(const TVTrackHit &,
                                   const TVector3   &xv) const
{
   return XvToMv(xv);
}

TVector3 EXITFBMeasLayer::HitToXv(const TVTrackHit &vht) const
{
   const EXITFBHit &ht = dynamic_cast<const EXITFBHit &>(vht);

   Double_t x = ht(0,0);
   Double_t y = ht(1,0);
   Double_t z = GetXc().Z();
   return TVector3(x, y, z);
}

void EXITFBMeasLayer::CalcDhDa(const TVTrackHit &vht,
                               const TVector3   &xxv,
                               const TKalMatrix &dxphiada,
                                     TKalMatrix &H) const
{
   // Calculate
   //    H = (@h/@a) = (@phi/@a, @z/@a)^t
   // where
   //        h(a) = (phi, z)^t: expected meas vector
   //        a = (drho, phi0, kappa, dz, tanl, t0)
   //

   Int_t sdim = H.GetNcols();
   Int_t hdim = 5;

   // Set H = (@h/@a) = (@phi/@a, @r/@a)^t
   
   for (Int_t i=0; i<hdim; i++) {
      H(0,i) = dxphiada(0,i);
      H(1,i) = dxphiada(1,i);
   }
   if (sdim == 6) {
      H(0,sdim-1) = 0.;
      H(1,sdim-1) = 0.;
   }
}

Int_t EXITFBMeasLayer::CalcXingPointWith(const TVTrack &hel,
                                               TVector3 &xx,
                                               Double_t &phi,
                                               Int_t     mode,
                                               Double_t  eps) const
{
   // This assumes nonzero B field.
   //
   // Copy helix parameters to local variables.
   //

   Double_t dr  = hel.GetDrho();        // drho
   Double_t fi0 = hel.GetPhi0();        // fi0
   Double_t dz  = hel.GetDz();          // dz 
   Double_t tnl = hel.GetTanLambda();   // tan lamda
   Double_t cpa = hel.GetKappa();       // kappa
   TVector3 X0  = hel.GetPivot();       // vector (x0, y0, z0)

   //
   // Check if charge is nonzero.
   //

   Int_t    chg = (Int_t)TMath::Sign(1.1,cpa);
   if (!chg) {
      cerr << ">>>> Error >>>> EXITFBMeasLayer::CalcXingPointWith" << endl
           << "      Kappa = 0 is invalid for a helix "            << endl;
      return -1;
   }

   //
   // Project everything to XY plane and calculate crossing points.
   //

   Double_t rho  = hel.GetRho();       // alpa/kappa
   Double_t csf0 = TMath::Cos(fi0);    // cos phi0
   Double_t snf0 = TMath::Sin(fi0);    // sin phi0
   Double_t z    = GetXc().Z();
  
   Double_t fi   = (-z + X0.Z() + dz) / (rho * tnl);
#if 1
   if (fi * chg * mode > 0. || 
       tnl * GetXc().Z() < 0. ||
       TMath::Abs(fi) >= 2*TMath::Pi()) return 0;
#endif
   Double_t x = X0.X() + dr*csf0 + rho*(csf0 - TMath::Cos(fi0 + fi));   // calculate X
   Double_t y = X0.Y() + dr*snf0 + rho*(snf0 - TMath::Sin(fi0 + fi));   // calculate Y
   xx.SetXYZ(x, y, z);

   if (IsOnSurface(xx)) {
      phi = fi;
      return 1;
   } else {
      return 0;
   }
}

void EXITFBMeasLayer::ProcessHit(const TVector3  &xx,
                                       TObjArray &hits)
{
   TKalMatrix h   = XvToMv(xx);
   Double_t   x = h(0, 0); 
   Double_t   y = h(1, 0); 

   Double_t dx = GetSigmaX();
   Double_t dy = GetSigmaY();
   x += gRandom->Gaus(0., dx);   // smearing x
   y += gRandom->Gaus(0., dy);   // smearing y

   Double_t meas [2];
   Double_t dmeas[2];
   meas [0] = x;
   meas [1] = y;
   dmeas[0] = dx;
   dmeas[1] = dy;
 
   Double_t b = EXITKalDetector::GetBfield();
   hits.Add(new EXITFBHit(*this, meas, dmeas, xx, b));
}
